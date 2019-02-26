% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:33
% EndTime: 2019-02-26 22:06:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
t109 = sin(pkin(6));
t111 = sin(qJ(2));
t124 = t109 * t111;
t112 = sin(qJ(1));
t123 = t109 * t112;
t113 = cos(qJ(2));
t122 = t109 * t113;
t114 = cos(qJ(1));
t121 = t109 * t114;
t120 = t112 * t111;
t119 = t112 * t113;
t118 = t114 * t111;
t117 = t114 * t113;
t110 = cos(pkin(6));
t103 = t110 * t118 + t119;
t108 = qJ(3) + pkin(11);
t106 = sin(t108);
t107 = cos(t108);
t116 = t103 * t106 + t107 * t121;
t115 = t103 * t107 - t106 * t121;
t105 = -t110 * t120 + t117;
t104 = t110 * t119 + t118;
t102 = t110 * t117 - t120;
t101 = t105 * t107 + t106 * t123;
t100 = t105 * t106 - t107 * t123;
t1 = [t102, t105, 0, 0, 0, 0; t104, t103, 0, 0, 0, 0; 0, t124, 0, 0, 0, 0; t115, t104 * t107, t100, 0, 0, 0; -t101, -t102 * t107, t116, 0, 0, 0; 0, -t107 * t122, t106 * t124 - t110 * t107, 0, 0, 0; -t116, -t104 * t106, t101, 0, 0, 0; t100, t102 * t106, t115, 0, 0, 0; 0, t106 * t122, t110 * t106 + t107 * t124, 0, 0, 0;];
JR_rot  = t1;
