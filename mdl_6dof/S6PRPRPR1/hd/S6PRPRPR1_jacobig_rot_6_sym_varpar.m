% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
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
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:28
% EndTime: 2019-02-26 19:46:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (24->11), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->18)
t122 = sin(pkin(10));
t123 = sin(pkin(6));
t132 = t122 * t123;
t125 = cos(pkin(10));
t131 = t125 * t123;
t121 = sin(pkin(11));
t124 = cos(pkin(11));
t127 = sin(qJ(2));
t128 = cos(qJ(2));
t130 = t128 * t121 + t127 * t124;
t129 = t127 * t121 - t128 * t124;
t126 = cos(pkin(6));
t120 = qJ(4) + pkin(12);
t119 = cos(t120);
t118 = sin(t120);
t115 = t130 * t126;
t114 = t129 * t126;
t1 = [0, t132, 0, -t122 * t114 + t125 * t130, 0 (-t122 * t115 - t125 * t129) * t118 - t119 * t132; 0, -t131, 0, t125 * t114 + t122 * t130, 0 (t125 * t115 - t122 * t129) * t118 + t119 * t131; 0, t126, 0, t129 * t123, 0, t130 * t118 * t123 - t126 * t119;];
Jg_rot  = t1;
