% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (118->26), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->24)
t102 = sin(qJ(6));
t105 = cos(qJ(6));
t101 = pkin(10) + qJ(4);
t100 = cos(t101);
t107 = cos(qJ(1));
t104 = sin(qJ(1));
t106 = cos(qJ(2));
t115 = t104 * t106;
t99 = sin(t101);
t93 = t107 * t100 + t99 * t115;
t94 = t100 * t115 - t107 * t99;
t113 = t94 * t102 - t93 * t105;
t103 = sin(qJ(2));
t110 = t100 * t105 + t102 * t99;
t117 = t110 * t103;
t114 = t107 * t106;
t112 = t102 * t93 + t105 * t94;
t95 = -t104 * t100 + t99 * t114;
t96 = t100 * t114 + t104 * t99;
t111 = t102 * t96 - t95 * t105;
t89 = t102 * t95 + t105 * t96;
t109 = t100 * t102 - t105 * t99;
t91 = t109 * t103;
t1 = [-t112, -t107 * t117, 0, t111, 0, -t111; t89, -t104 * t117, 0, t113, 0, -t113; 0, t110 * t106, 0, t91, 0, -t91; t113, t107 * t91, 0, t89, 0, -t89; -t111, t104 * t91, 0, t112, 0, -t112; 0, -t109 * t106, 0, t117, 0, -t117; t104 * t103, -t114, 0, 0, 0, 0; -t107 * t103, -t115, 0, 0, 0, 0; 0, -t103, 0, 0, 0, 0;];
JR_rot  = t1;
