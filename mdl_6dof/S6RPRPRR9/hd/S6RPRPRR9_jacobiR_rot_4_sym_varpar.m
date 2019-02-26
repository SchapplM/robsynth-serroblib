% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR9_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (70->22), mult. (202->42), div. (0->0), fcn. (284->12), ass. (0->33)
t119 = sin(pkin(6));
t125 = sin(qJ(1));
t137 = t119 * t125;
t127 = cos(qJ(1));
t136 = t119 * t127;
t117 = sin(pkin(12));
t135 = t125 * t117;
t121 = cos(pkin(12));
t134 = t125 * t121;
t133 = t127 * t117;
t132 = t127 * t121;
t116 = sin(pkin(13));
t120 = cos(pkin(13));
t124 = sin(qJ(3));
t126 = cos(qJ(3));
t131 = t126 * t116 + t124 * t120;
t112 = t124 * t116 - t126 * t120;
t118 = sin(pkin(7));
t104 = t112 * t118;
t122 = cos(pkin(7));
t106 = t112 * t122;
t123 = cos(pkin(6));
t108 = -t123 * t132 + t135;
t109 = t123 * t133 + t134;
t130 = -t104 * t136 - t108 * t106 + t109 * t131;
t105 = t131 * t118;
t107 = t131 * t122;
t129 = t105 * t136 + t108 * t107 + t109 * t112;
t110 = -t123 * t134 - t133;
t111 = -t123 * t135 + t132;
t128 = t105 * t137 + t110 * t107 - t111 * t112;
t103 = -t104 * t137 - t110 * t106 - t111 * t131;
t1 = [t129, 0, t103, 0, 0, 0; t128, 0, -t130, 0, 0, 0; 0, 0, -t123 * t104 + (-t106 * t121 - t117 * t131) * t119, 0, 0, 0; t130, 0, -t128, 0, 0, 0; t103, 0, t129, 0, 0, 0; 0, 0, -t123 * t105 + (-t107 * t121 + t112 * t117) * t119, 0, 0, 0; -t108 * t118 + t122 * t136, 0, 0, 0, 0, 0; -t110 * t118 + t122 * t137, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
