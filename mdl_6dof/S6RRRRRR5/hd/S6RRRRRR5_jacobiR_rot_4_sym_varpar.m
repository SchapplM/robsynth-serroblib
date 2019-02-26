% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR5_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:49:26
% EndTime: 2019-02-26 22:49:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (77->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
t108 = sin(pkin(6));
t110 = sin(qJ(2));
t122 = t108 * t110;
t111 = sin(qJ(1));
t121 = t108 * t111;
t112 = cos(qJ(2));
t120 = t108 * t112;
t113 = cos(qJ(1));
t119 = t108 * t113;
t118 = t111 * t110;
t117 = t111 * t112;
t116 = t113 * t110;
t115 = t113 * t112;
t109 = cos(pkin(6));
t101 = t109 * t116 + t117;
t107 = qJ(3) + qJ(4);
t105 = sin(t107);
t106 = cos(t107);
t95 = -t101 * t106 + t105 * t119;
t114 = t101 * t105 + t106 * t119;
t103 = -t109 * t118 + t115;
t102 = t109 * t117 + t116;
t100 = t109 * t115 - t118;
t99 = -t109 * t105 - t106 * t122;
t98 = -t105 * t122 + t109 * t106;
t97 = t103 * t106 + t105 * t121;
t96 = -t103 * t105 + t106 * t121;
t1 = [t95, -t102 * t106, t96, t96, 0, 0; t97, t100 * t106, -t114, -t114, 0, 0; 0, t106 * t120, t98, t98, 0, 0; t114, t102 * t105, -t97, -t97, 0, 0; t96, -t100 * t105, t95, t95, 0, 0; 0, -t105 * t120, t99, t99, 0, 0; t100, t103, 0, 0, 0, 0; t102, t101, 0, 0, 0, 0; 0, t122, 0, 0, 0, 0;];
JR_rot  = t1;
