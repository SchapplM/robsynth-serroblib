% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:16
% EndTime: 2019-02-26 19:48:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->20), mult. (107->52), div. (0->0), fcn. (156->10), ass. (0->27)
t105 = pkin(11) + qJ(4);
t104 = cos(t105);
t106 = sin(pkin(12));
t121 = t104 * t106;
t109 = cos(pkin(12));
t120 = t104 * t109;
t113 = cos(qJ(2));
t119 = t104 * t113;
t107 = sin(pkin(10));
t108 = sin(pkin(6));
t118 = t107 * t108;
t110 = cos(pkin(10));
t117 = t108 * t110;
t112 = sin(qJ(2));
t116 = t108 * t112;
t111 = cos(pkin(6));
t115 = t111 * t112;
t114 = t111 * t113;
t103 = sin(t105);
t102 = -t107 * t115 + t110 * t113;
t101 = -t107 * t114 - t110 * t112;
t100 = t107 * t113 + t110 * t115;
t99 = -t107 * t112 + t110 * t114;
t98 = -t103 * t116 + t111 * t104;
t97 = -t102 * t103 + t104 * t118;
t96 = -t100 * t103 - t104 * t117;
t1 = [0, t101 * t120 + t102 * t106, 0, t97 * t109, 0, 0; 0, t100 * t106 + t99 * t120, 0, t96 * t109, 0, 0; 0 (t106 * t112 + t109 * t119) * t108, 0, t98 * t109, 0, 0; 0, -t101 * t121 + t102 * t109, 0, -t97 * t106, 0, 0; 0, t100 * t109 - t99 * t121, 0, -t96 * t106, 0, 0; 0 (-t106 * t119 + t109 * t112) * t108, 0, -t98 * t106, 0, 0; 0, t101 * t103, 0, t102 * t104 + t103 * t118, 0, 0; 0, t99 * t103, 0, t100 * t104 - t103 * t117, 0, 0; 0, t108 * t113 * t103, 0, t111 * t103 + t104 * t116, 0, 0;];
JR_rot  = t1;
