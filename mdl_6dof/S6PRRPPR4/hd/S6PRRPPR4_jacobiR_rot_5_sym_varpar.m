% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:44
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:44:03
% EndTime: 2019-02-22 09:44:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (31->16), mult. (107->50), div. (0->0), fcn. (156->10), ass. (0->25)
t111 = sin(pkin(11));
t119 = cos(qJ(3));
t127 = t111 * t119;
t113 = sin(pkin(6));
t117 = sin(qJ(3));
t126 = t113 * t117;
t125 = t113 * t119;
t114 = cos(pkin(11));
t124 = t114 * t119;
t116 = cos(pkin(6));
t118 = sin(qJ(2));
t123 = t116 * t118;
t120 = cos(qJ(2));
t122 = t116 * t120;
t121 = t119 * t120;
t115 = cos(pkin(10));
t112 = sin(pkin(10));
t110 = t116 * t119 - t118 * t126;
t109 = -t112 * t123 + t115 * t120;
t108 = -t112 * t122 - t115 * t118;
t107 = t112 * t120 + t115 * t123;
t106 = -t112 * t118 + t115 * t122;
t105 = -t109 * t117 + t112 * t125;
t104 = -t107 * t117 - t115 * t125;
t1 = [0, t108 * t124 + t109 * t111, t105 * t114, 0, 0, 0; 0, t106 * t124 + t107 * t111, t104 * t114, 0, 0, 0; 0 (t111 * t118 + t114 * t121) * t113, t110 * t114, 0, 0, 0; 0, t108 * t117, t109 * t119 + t112 * t126, 0, 0, 0; 0, t106 * t117, t107 * t119 - t115 * t126, 0, 0, 0; 0, t120 * t126, t116 * t117 + t118 * t125, 0, 0, 0; 0, t108 * t127 - t109 * t114, t105 * t111, 0, 0, 0; 0, t106 * t127 - t107 * t114, t104 * t111, 0, 0, 0; 0 (t111 * t121 - t114 * t118) * t113, t110 * t111, 0, 0, 0;];
JR_rot  = t1;
