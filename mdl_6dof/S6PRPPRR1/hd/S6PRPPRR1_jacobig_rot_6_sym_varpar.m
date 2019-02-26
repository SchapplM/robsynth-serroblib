% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPPRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->11), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->18)
t119 = sin(pkin(10));
t120 = sin(pkin(6));
t129 = t119 * t120;
t122 = cos(pkin(10));
t128 = t122 * t120;
t118 = sin(pkin(11));
t121 = cos(pkin(11));
t124 = sin(qJ(2));
t125 = cos(qJ(2));
t127 = t125 * t118 + t124 * t121;
t126 = t124 * t118 - t125 * t121;
t123 = cos(pkin(6));
t117 = pkin(12) + qJ(5);
t116 = cos(t117);
t115 = sin(t117);
t112 = t127 * t123;
t111 = t126 * t123;
t1 = [0, t129, 0, 0, -t119 * t111 + t122 * t127 (-t119 * t112 - t122 * t126) * t115 - t116 * t129; 0, -t128, 0, 0, t122 * t111 + t119 * t127 (t122 * t112 - t119 * t126) * t115 + t116 * t128; 0, t123, 0, 0, t126 * t120, t127 * t115 * t120 - t123 * t116;];
Jg_rot  = t1;
