% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:04
% EndTime: 2019-02-26 21:06:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (48->22), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->27)
t103 = cos(qJ(6));
t100 = sin(qJ(4));
t101 = sin(qJ(3));
t106 = cos(qJ(1));
t113 = t106 * t101;
t102 = sin(qJ(1));
t104 = cos(qJ(4));
t115 = t102 * t104;
t95 = t100 * t113 + t115;
t112 = t106 * t104;
t96 = -t102 * t100 + t101 * t112;
t99 = sin(qJ(6));
t110 = t95 * t103 - t96 * t99;
t116 = t102 * t101;
t105 = cos(qJ(3));
t114 = t102 * t105;
t111 = t106 * t105;
t93 = t100 * t116 - t112;
t94 = t106 * t100 + t101 * t115;
t89 = t94 * t103 + t93 * t99;
t88 = t93 * t103 - t94 * t99;
t109 = t96 * t103 + t95 * t99;
t108 = t100 * t99 + t103 * t104;
t107 = t100 * t103 - t104 * t99;
t92 = t108 * t105;
t91 = t107 * t105;
t1 = [t109, 0, t102 * t92, -t88, 0, t88; t89, 0, -t108 * t111, t110, 0, -t110; 0, 0, -t108 * t101, -t91, 0, t91; t110, 0, t107 * t114, t89, 0, -t89; t88, 0, -t106 * t91, -t109, 0, t109; 0, 0, -t107 * t101, t92, 0, -t92; t111, 0, -t116, 0, 0, 0; t114, 0, t113, 0, 0, 0; 0, 0, -t105, 0, 0, 0;];
JR_rot  = t1;
