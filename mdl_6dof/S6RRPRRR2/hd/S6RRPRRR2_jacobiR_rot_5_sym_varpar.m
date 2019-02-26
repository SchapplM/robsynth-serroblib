% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR2
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

function JR_rot = S6RRPRRR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:40
% EndTime: 2019-02-26 21:54:40
% DurationCPUTime: 0.04s
% Computational Cost: add. (77->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t103 = qJ(2) + pkin(11) + qJ(4);
t102 = cos(t103);
t104 = sin(qJ(5));
t114 = t102 * t104;
t105 = sin(qJ(1));
t113 = t105 * t104;
t106 = cos(qJ(5));
t112 = t105 * t106;
t107 = cos(qJ(1));
t111 = t107 * t104;
t110 = t107 * t106;
t101 = sin(t103);
t109 = t101 * t112;
t108 = t101 * t110;
t100 = t107 * t102;
t99 = t102 * t106;
t98 = t105 * t102;
t97 = t101 * t111;
t96 = t101 * t113;
t95 = t102 * t110 + t113;
t94 = -t102 * t111 + t112;
t93 = -t102 * t112 + t111;
t92 = t102 * t113 + t110;
t1 = [t93, -t108, 0, -t108, t94, 0; t95, -t109, 0, -t109, -t92, 0; 0, t99, 0, t99, -t101 * t104, 0; t92, t97, 0, t97, -t95, 0; t94, t96, 0, t96, t93, 0; 0, -t114, 0, -t114, -t101 * t106, 0; -t105 * t101, t100, 0, t100, 0, 0; t107 * t101, t98, 0, t98, 0, 0; 0, t101, 0, t101, 0, 0;];
JR_rot  = t1;
