% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (80->16), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t115 = pkin(10) + qJ(3);
t111 = sin(t115);
t116 = qJ(4) + qJ(5);
t113 = sin(t116);
t123 = t111 * t113;
t117 = sin(qJ(1));
t122 = t117 * t113;
t114 = cos(t116);
t121 = t117 * t114;
t118 = cos(qJ(1));
t120 = t118 * t113;
t119 = t118 * t114;
t112 = cos(t115);
t109 = t111 * t114;
t108 = t112 * t119 + t122;
t107 = t112 * t120 - t121;
t106 = t112 * t121 - t120;
t105 = -t112 * t122 - t119;
t1 = [-t106, 0, -t111 * t119, -t107, -t107, 0; t108, 0, -t111 * t121, t105, t105, 0; 0, 0, t112 * t114, -t123, -t123, 0; -t117 * t111, 0, t118 * t112, 0, 0, 0; t118 * t111, 0, t117 * t112, 0, 0, 0; 0, 0, t111, 0, 0, 0; t105, 0, -t111 * t120, t108, t108, 0; t107, 0, -t111 * t122, t106, t106, 0; 0, 0, t112 * t113, t109, t109, 0;];
JR_rot  = t1;
