% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (137->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t110 = qJ(2) + qJ(3) + pkin(11) + qJ(5);
t109 = cos(t110);
t111 = sin(qJ(6));
t121 = t109 * t111;
t112 = sin(qJ(1));
t120 = t112 * t111;
t113 = cos(qJ(6));
t119 = t112 * t113;
t114 = cos(qJ(1));
t118 = t114 * t111;
t117 = t114 * t113;
t108 = sin(t110);
t116 = t108 * t119;
t115 = t108 * t117;
t107 = t114 * t109;
t106 = t109 * t113;
t105 = t112 * t109;
t104 = t108 * t118;
t103 = t108 * t120;
t102 = t109 * t117 + t120;
t101 = -t109 * t118 + t119;
t100 = -t109 * t119 + t118;
t99 = t109 * t120 + t117;
t1 = [t100, -t115, -t115, 0, -t115, t101; t102, -t116, -t116, 0, -t116, -t99; 0, t106, t106, 0, t106, -t108 * t111; t99, t104, t104, 0, t104, -t102; t101, t103, t103, 0, t103, t100; 0, -t121, -t121, 0, -t121, -t108 * t113; -t112 * t108, t107, t107, 0, t107, 0; t114 * t108, t105, t105, 0, t105, 0; 0, t108, t108, 0, t108, 0;];
JR_rot  = t1;
