% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:13
% EndTime: 2019-02-26 20:49:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
t108 = qJ(3) + pkin(11) + qJ(5);
t105 = cos(t108);
t110 = sin(qJ(6));
t118 = t105 * t110;
t109 = qJ(1) + pkin(10);
t106 = sin(t109);
t117 = t106 * t110;
t111 = cos(qJ(6));
t116 = t106 * t111;
t107 = cos(t109);
t115 = t107 * t110;
t114 = t107 * t111;
t104 = sin(t108);
t113 = t104 * t116;
t112 = t104 * t114;
t103 = t105 * t111;
t102 = t107 * t105;
t101 = t106 * t105;
t100 = t104 * t115;
t99 = t104 * t117;
t98 = t105 * t114 + t117;
t97 = -t105 * t115 + t116;
t96 = -t105 * t116 + t115;
t95 = t105 * t117 + t114;
t1 = [t96, 0, -t112, 0, -t112, t97; t98, 0, -t113, 0, -t113, -t95; 0, 0, t103, 0, t103, -t104 * t110; t95, 0, t100, 0, t100, -t98; t97, 0, t99, 0, t99, t96; 0, 0, -t118, 0, -t118, -t104 * t111; -t106 * t104, 0, t102, 0, t102, 0; t107 * t104, 0, t101, 0, t101, 0; 0, 0, t104, 0, t104, 0;];
JR_rot  = t1;
