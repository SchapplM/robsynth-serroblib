% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:26
% EndTime: 2019-02-26 21:16:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (137->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t107 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
t106 = cos(t107);
t108 = sin(qJ(6));
t118 = t106 * t108;
t109 = sin(qJ(1));
t117 = t109 * t108;
t110 = cos(qJ(6));
t116 = t109 * t110;
t111 = cos(qJ(1));
t115 = t111 * t108;
t114 = t111 * t110;
t105 = sin(t107);
t113 = t105 * t116;
t112 = t105 * t114;
t104 = t111 * t106;
t103 = t106 * t110;
t102 = t109 * t106;
t101 = t105 * t115;
t100 = t105 * t117;
t99 = t106 * t114 + t117;
t98 = -t106 * t115 + t116;
t97 = -t106 * t116 + t115;
t96 = t106 * t117 + t114;
t1 = [t97, 0, -t112, -t112, -t112, t98; t99, 0, -t113, -t113, -t113, -t96; 0, 0, t103, t103, t103, -t105 * t108; t96, 0, t101, t101, t101, -t99; t98, 0, t100, t100, t100, t97; 0, 0, -t118, -t118, -t118, -t105 * t110; -t109 * t105, 0, t104, t104, t104, 0; t111 * t105, 0, t102, t102, t102, 0; 0, 0, t105, t105, t105, 0;];
JR_rot  = t1;
