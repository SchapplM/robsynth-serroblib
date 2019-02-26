% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPPRRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:54
% EndTime: 2019-02-26 19:38:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (62->25), mult. (184->54), div. (0->0), fcn. (252->14), ass. (0->34)
t95 = sin(pkin(13));
t99 = sin(pkin(6));
t119 = t95 * t99;
t96 = sin(pkin(12));
t118 = t96 * t95;
t98 = sin(pkin(7));
t117 = t98 * t99;
t104 = cos(pkin(7));
t116 = t104 * t99;
t101 = cos(pkin(13));
t115 = t96 * t101;
t102 = cos(pkin(12));
t105 = cos(pkin(6));
t114 = t102 * t105;
t100 = cos(pkin(14));
t103 = cos(pkin(8));
t90 = t101 * t114 - t118;
t109 = -t102 * t117 + t104 * t90;
t91 = t95 * t114 + t115;
t94 = sin(pkin(14));
t97 = sin(pkin(8));
t113 = t103 * (t109 * t100 - t91 * t94) + (-t102 * t116 - t90 * t98) * t97;
t92 = -t102 * t95 - t105 * t115;
t110 = t104 * t92 + t96 * t117;
t93 = t102 * t101 - t105 * t118;
t112 = t103 * (t110 * t100 - t93 * t94) + (t96 * t116 - t92 * t98) * t97;
t108 = t101 * t116 + t105 * t98;
t111 = t103 * (t108 * t100 - t94 * t119) + (-t101 * t117 + t105 * t104) * t97;
t107 = cos(qJ(4));
t106 = sin(qJ(4));
t86 = t100 * t119 + t108 * t94;
t84 = t93 * t100 + t110 * t94;
t82 = t91 * t100 + t109 * t94;
t1 = [0, 0, 0, -t84 * t106 + t112 * t107, 0, 0; 0, 0, 0, -t82 * t106 + t113 * t107, 0, 0; 0, 0, 0, -t86 * t106 + t111 * t107, 0, 0; 0, 0, 0, -t112 * t106 - t84 * t107, 0, 0; 0, 0, 0, -t113 * t106 - t82 * t107, 0, 0; 0, 0, 0, -t111 * t106 - t86 * t107, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
