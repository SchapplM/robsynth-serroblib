% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:18
% EndTime: 2019-02-26 19:50:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t122 = sin(pkin(10));
t123 = sin(pkin(6));
t134 = t122 * t123;
t125 = cos(pkin(10));
t133 = t125 * t123;
t121 = sin(pkin(11));
t124 = cos(pkin(11));
t128 = sin(qJ(2));
t130 = cos(qJ(2));
t132 = t130 * t121 + t128 * t124;
t131 = t128 * t121 - t130 * t124;
t129 = cos(qJ(4));
t127 = sin(qJ(4));
t126 = cos(pkin(6));
t118 = t132 * t126;
t117 = t131 * t126;
t1 = [0, t134, 0, -t122 * t117 + t125 * t132 (-t122 * t118 - t125 * t131) * t127 - t129 * t134, 0; 0, -t133, 0, t125 * t117 + t122 * t132 (t125 * t118 - t122 * t131) * t127 + t129 * t133, 0; 0, t126, 0, t131 * t123, t132 * t127 * t123 - t126 * t129, 0;];
Jg_rot  = t1;
