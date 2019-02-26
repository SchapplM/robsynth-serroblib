% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:43
% EndTime: 2019-02-26 20:11:43
% DurationCPUTime: 0.05s
% Computational Cost: add. (59->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
t116 = sin(pkin(11));
t117 = sin(pkin(6));
t127 = t116 * t117;
t118 = cos(pkin(11));
t126 = t117 * t118;
t120 = sin(qJ(2));
t125 = t117 * t120;
t121 = cos(qJ(2));
t124 = t117 * t121;
t119 = cos(pkin(6));
t123 = t119 * t120;
t122 = t119 * t121;
t115 = qJ(3) + qJ(4);
t114 = cos(t115);
t113 = sin(t115);
t112 = -t116 * t123 + t118 * t121;
t111 = -t116 * t122 - t118 * t120;
t110 = t116 * t121 + t118 * t123;
t109 = -t116 * t120 + t118 * t122;
t108 = t119 * t113 + t114 * t125;
t107 = t113 * t125 - t119 * t114;
t106 = t112 * t114 + t113 * t127;
t105 = t112 * t113 - t114 * t127;
t104 = t110 * t114 - t113 * t126;
t103 = t110 * t113 + t114 * t126;
t1 = [0, t112, 0, 0, 0, 0; 0, t110, 0, 0, 0, 0; 0, t125, 0, 0, 0, 0; 0, -t111 * t114, t105, t105, 0, 0; 0, -t109 * t114, t103, t103, 0, 0; 0, -t114 * t124, t107, t107, 0, 0; 0, t111 * t113, t106, t106, 0, 0; 0, t109 * t113, t104, t104, 0, 0; 0, t113 * t124, t108, t108, 0, 0;];
JR_rot  = t1;
