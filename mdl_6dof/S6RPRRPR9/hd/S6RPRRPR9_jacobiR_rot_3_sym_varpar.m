% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR9_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiR_rot_3_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:28
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (40->17), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->27)
t103 = cos(qJ(1));
t96 = sin(pkin(6));
t111 = t103 * t96;
t101 = sin(qJ(1));
t99 = cos(pkin(6));
t110 = t103 * t99;
t94 = sin(pkin(12));
t97 = cos(pkin(12));
t88 = t101 * t94 - t97 * t110;
t95 = sin(pkin(7));
t98 = cos(pkin(7));
t117 = t95 * t111 + t88 * t98;
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t89 = t101 * t97 + t94 * t110;
t116 = t89 * t100 + t117 * t102;
t114 = t94 * t96;
t113 = t101 * t96;
t112 = t101 * t99;
t107 = t96 * t97 * t98 + t95 * t99;
t90 = -t103 * t94 - t97 * t112;
t106 = t95 * t113 + t90 * t98;
t104 = t117 * t100 - t89 * t102;
t91 = t103 * t97 - t94 * t112;
t87 = t106 * t100 + t91 * t102;
t86 = -t91 * t100 + t106 * t102;
t1 = [t104, 0, t86, 0, 0, 0; t87, 0, -t116, 0, 0, 0; 0, 0, -t100 * t114 + t107 * t102, 0, 0, 0; t116, 0, -t87, 0, 0, 0; t86, 0, t104, 0, 0, 0; 0, 0, -t107 * t100 - t102 * t114, 0, 0, 0; t98 * t111 - t88 * t95, 0, 0, 0, 0, 0; t98 * t113 - t90 * t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
