% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPPRR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->14), mult. (50->32), div. (0->0), fcn. (76->10), ass. (0->19)
t115 = sin(pkin(10));
t116 = sin(pkin(6));
t127 = t115 * t116;
t118 = cos(pkin(10));
t126 = t118 * t116;
t119 = cos(pkin(6));
t121 = sin(qJ(2));
t125 = t119 * t121;
t123 = cos(qJ(2));
t124 = t119 * t123;
t122 = cos(qJ(5));
t120 = sin(qJ(5));
t117 = cos(pkin(11));
t114 = sin(pkin(11));
t113 = -t115 * t125 + t118 * t123;
t112 = t115 * t124 + t118 * t121;
t111 = t115 * t123 + t118 * t125;
t110 = t115 * t121 - t118 * t124;
t1 = [0, t127, 0, 0, -t112 * t117 + t113 * t114 (t112 * t114 + t113 * t117) * t120 + t122 * t127; 0, -t126, 0, 0, -t110 * t117 + t111 * t114 (t110 * t114 + t111 * t117) * t120 - t122 * t126; 0, t119, 0, 0 (t114 * t121 + t117 * t123) * t116, t119 * t122 + (-t114 * t123 + t117 * t121) * t120 * t116;];
Jg_rot  = t1;
