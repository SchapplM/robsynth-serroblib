% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:36
% EndTime: 2019-02-26 22:10:36
% DurationCPUTime: 0.04s
% Computational Cost: add. (78->18), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
t115 = qJ(2) + qJ(3);
t112 = sin(t115);
t116 = sin(qJ(1));
t123 = t116 * t112;
t113 = cos(t115);
t108 = t116 * t113;
t117 = cos(qJ(1));
t122 = t117 * t112;
t109 = t117 * t113;
t114 = pkin(10) + qJ(5);
t110 = sin(t114);
t121 = t110 * t123;
t111 = cos(t114);
t120 = t111 * t123;
t119 = t110 * t122;
t118 = t111 * t122;
t107 = t113 * t111;
t106 = t113 * t110;
t105 = t111 * t109 + t116 * t110;
t104 = t110 * t109 - t116 * t111;
t103 = t111 * t108 - t117 * t110;
t102 = -t110 * t108 - t117 * t111;
t1 = [-t103, -t118, -t118, 0, -t104, 0; t105, -t120, -t120, 0, t102, 0; 0, t107, t107, 0, -t112 * t110, 0; -t123, t109, t109, 0, 0, 0; t122, t108, t108, 0, 0, 0; 0, t112, t112, 0, 0, 0; t102, -t119, -t119, 0, t105, 0; t104, -t121, -t121, 0, t103, 0; 0, t106, t106, 0, t112 * t111, 0;];
JR_rot  = t1;
