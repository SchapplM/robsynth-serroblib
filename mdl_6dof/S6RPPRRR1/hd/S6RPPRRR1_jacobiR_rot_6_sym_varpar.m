% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:43
% EndTime: 2019-02-26 20:34:43
% DurationCPUTime: 0.04s
% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
t106 = pkin(11) + qJ(4) + qJ(5);
t103 = cos(t106);
t108 = sin(qJ(6));
t116 = t103 * t108;
t107 = qJ(1) + pkin(10);
t104 = sin(t107);
t115 = t104 * t108;
t109 = cos(qJ(6));
t114 = t104 * t109;
t105 = cos(t107);
t113 = t105 * t108;
t112 = t105 * t109;
t102 = sin(t106);
t111 = t102 * t114;
t110 = t102 * t112;
t101 = t103 * t109;
t100 = t105 * t103;
t99 = t104 * t103;
t98 = t102 * t113;
t97 = t102 * t115;
t96 = t103 * t112 + t115;
t95 = -t103 * t113 + t114;
t94 = -t103 * t114 + t113;
t93 = t103 * t115 + t112;
t1 = [t94, 0, 0, -t110, -t110, t95; t96, 0, 0, -t111, -t111, -t93; 0, 0, 0, t101, t101, -t102 * t108; t93, 0, 0, t98, t98, -t96; t95, 0, 0, t97, t97, t94; 0, 0, 0, -t116, -t116, -t102 * t109; -t104 * t102, 0, 0, t100, t100, 0; t105 * t102, 0, 0, t99, t99, 0; 0, 0, 0, t102, t102, 0;];
JR_rot  = t1;
