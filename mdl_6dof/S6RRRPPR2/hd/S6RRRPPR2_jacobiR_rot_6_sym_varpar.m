% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:56
% EndTime: 2019-02-26 22:03:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (74->13), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t109 = qJ(2) + qJ(3) + pkin(10);
t107 = sin(t109);
t111 = sin(qJ(1));
t119 = t111 * t107;
t110 = sin(qJ(6));
t118 = t111 * t110;
t112 = cos(qJ(6));
t117 = t111 * t112;
t113 = cos(qJ(1));
t116 = t113 * t107;
t115 = t113 * t110;
t114 = t113 * t112;
t108 = cos(t109);
t106 = t107 * t112;
t105 = t107 * t110;
t104 = t108 * t114;
t103 = t108 * t115;
t102 = t108 * t117;
t101 = t108 * t118;
t100 = -t107 * t118 + t114;
t99 = t107 * t117 + t115;
t98 = t107 * t115 + t117;
t97 = t107 * t114 - t118;
t1 = [t100, t103, t103, 0, 0, t97; t98, t101, t101, 0, 0, t99; 0, t105, t105, 0, 0, -t108 * t112; -t99, t104, t104, 0, 0, -t98; t97, t102, t102, 0, 0, t100; 0, t106, t106, 0, 0, t108 * t110; -t111 * t108, -t116, -t116, 0, 0, 0; t113 * t108, -t119, -t119, 0, 0, 0; 0, t108, t108, 0, 0, 0;];
JR_rot  = t1;
