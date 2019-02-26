% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:06
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.05s
% Computational Cost: add. (89->14), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->28)
t96 = sin(pkin(11));
t97 = sin(pkin(6));
t109 = t96 * t97;
t98 = cos(pkin(11));
t108 = t97 * t98;
t101 = cos(qJ(2));
t107 = t101 * t97;
t100 = sin(qJ(2));
t106 = t96 * t100;
t105 = t96 * t101;
t104 = t97 * t100;
t103 = t98 * t100;
t102 = t98 * t101;
t99 = cos(pkin(6));
t95 = qJ(3) + pkin(12) + qJ(5);
t94 = cos(t95);
t93 = sin(t95);
t92 = -t99 * t106 + t102;
t91 = -t99 * t105 - t103;
t90 = t99 * t103 + t105;
t89 = t99 * t102 - t106;
t88 = -t94 * t104 - t99 * t93;
t87 = -t93 * t104 + t99 * t94;
t86 = -t93 * t109 - t92 * t94;
t85 = t94 * t109 - t92 * t93;
t84 = t93 * t108 - t90 * t94;
t83 = -t94 * t108 - t90 * t93;
t1 = [0, t91 * t94, t85, 0, t85, 0; 0, t89 * t94, t83, 0, t83, 0; 0, t94 * t107, t87, 0, t87, 0; 0, -t91 * t93, t86, 0, t86, 0; 0, -t89 * t93, t84, 0, t84, 0; 0, -t93 * t107, t88, 0, t88, 0; 0, t92, 0, 0, 0, 0; 0, t90, 0, 0, 0, 0; 0, t104, 0, 0, 0, 0;];
JR_rot  = t1;
