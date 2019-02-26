% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t119 = sin(pkin(10));
t120 = sin(pkin(6));
t131 = t119 * t120;
t122 = cos(pkin(10));
t130 = t122 * t120;
t118 = sin(pkin(11));
t121 = cos(pkin(11));
t125 = sin(qJ(2));
t127 = cos(qJ(2));
t129 = t127 * t118 + t125 * t121;
t128 = t125 * t118 - t127 * t121;
t126 = cos(qJ(4));
t124 = sin(qJ(4));
t123 = cos(pkin(6));
t115 = t129 * t123;
t114 = t128 * t123;
t1 = [0, t131, 0, -t119 * t114 + t122 * t129, 0 (-t119 * t115 - t122 * t128) * t124 - t126 * t131; 0, -t130, 0, t122 * t114 + t119 * t129, 0 (t122 * t115 - t119 * t128) * t124 + t126 * t130; 0, t123, 0, t128 * t120, 0, t129 * t124 * t120 - t123 * t126;];
Jg_rot  = t1;
