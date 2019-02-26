% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:39
% EndTime: 2019-02-26 19:47:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t114 = sin(pkin(10));
t115 = sin(pkin(6));
t126 = t114 * t115;
t117 = cos(pkin(10));
t125 = t117 * t115;
t113 = sin(pkin(11));
t116 = cos(pkin(11));
t120 = sin(qJ(2));
t122 = cos(qJ(2));
t124 = t122 * t113 + t120 * t116;
t123 = t120 * t113 - t122 * t116;
t121 = cos(qJ(4));
t119 = sin(qJ(4));
t118 = cos(pkin(6));
t110 = t124 * t118;
t109 = t123 * t118;
t1 = [0, t126, 0, -t114 * t109 + t117 * t124, 0 (-t114 * t110 - t117 * t123) * t121 + t119 * t126; 0, -t125, 0, t117 * t109 + t114 * t124, 0 (t117 * t110 - t114 * t123) * t121 - t119 * t125; 0, t118, 0, t123 * t115, 0, t124 * t121 * t115 + t118 * t119;];
Jg_rot  = t1;
