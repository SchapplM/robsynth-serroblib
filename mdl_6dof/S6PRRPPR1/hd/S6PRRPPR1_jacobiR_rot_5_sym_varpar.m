% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:42
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:42:14
% EndTime: 2019-02-22 09:42:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (61->20), mult. (107->52), div. (0->0), fcn. (156->10), ass. (0->27)
t108 = qJ(3) + pkin(11);
t107 = cos(t108);
t109 = sin(pkin(12));
t124 = t107 * t109;
t112 = cos(pkin(12));
t123 = t107 * t112;
t116 = cos(qJ(2));
t122 = t107 * t116;
t110 = sin(pkin(10));
t111 = sin(pkin(6));
t121 = t110 * t111;
t113 = cos(pkin(10));
t120 = t111 * t113;
t115 = sin(qJ(2));
t119 = t111 * t115;
t114 = cos(pkin(6));
t118 = t114 * t115;
t117 = t114 * t116;
t106 = sin(t108);
t105 = -t110 * t118 + t113 * t116;
t104 = -t110 * t117 - t113 * t115;
t103 = t110 * t116 + t113 * t118;
t102 = -t110 * t115 + t113 * t117;
t101 = -t106 * t119 + t114 * t107;
t100 = -t105 * t106 + t107 * t121;
t99 = -t103 * t106 - t107 * t120;
t1 = [0, t104 * t123 + t105 * t109, t100 * t112, 0, 0, 0; 0, t102 * t123 + t103 * t109, t99 * t112, 0, 0, 0; 0 (t109 * t115 + t112 * t122) * t111, t101 * t112, 0, 0, 0; 0, -t104 * t124 + t105 * t112, -t100 * t109, 0, 0, 0; 0, -t102 * t124 + t103 * t112, -t99 * t109, 0, 0, 0; 0 (-t109 * t122 + t112 * t115) * t111, -t101 * t109, 0, 0, 0; 0, t104 * t106, t105 * t107 + t106 * t121, 0, 0, 0; 0, t102 * t106, t103 * t107 - t106 * t120, 0, 0, 0; 0, t111 * t116 * t106, t114 * t106 + t107 * t119, 0, 0, 0;];
JR_rot  = t1;
