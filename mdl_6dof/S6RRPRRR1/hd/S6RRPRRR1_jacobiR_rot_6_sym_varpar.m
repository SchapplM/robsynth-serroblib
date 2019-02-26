% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:01
% EndTime: 2019-02-26 21:54:01
% DurationCPUTime: 0.04s
% Computational Cost: add. (137->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t109 = qJ(2) + pkin(11) + qJ(4) + qJ(5);
t108 = cos(t109);
t110 = sin(qJ(6));
t120 = t108 * t110;
t111 = sin(qJ(1));
t119 = t111 * t110;
t112 = cos(qJ(6));
t118 = t111 * t112;
t113 = cos(qJ(1));
t117 = t113 * t110;
t116 = t113 * t112;
t107 = sin(t109);
t115 = t107 * t118;
t114 = t107 * t116;
t106 = t113 * t108;
t105 = t108 * t112;
t104 = t111 * t108;
t103 = t107 * t117;
t102 = t107 * t119;
t101 = t108 * t116 + t119;
t100 = -t108 * t117 + t118;
t99 = -t108 * t118 + t117;
t98 = t108 * t119 + t116;
t1 = [t99, -t114, 0, -t114, -t114, t100; t101, -t115, 0, -t115, -t115, -t98; 0, t105, 0, t105, t105, -t107 * t110; t98, t103, 0, t103, t103, -t101; t100, t102, 0, t102, t102, t99; 0, -t120, 0, -t120, -t120, -t107 * t112; -t111 * t107, t106, 0, t106, t106, 0; t113 * t107, t104, 0, t104, t104, 0; 0, t107, 0, t107, t107, 0;];
JR_rot  = t1;
