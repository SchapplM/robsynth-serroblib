% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP2
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
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP2_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:55
% EndTime: 2019-02-26 19:50:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t111 = sin(pkin(10));
t112 = sin(pkin(6));
t123 = t111 * t112;
t114 = cos(pkin(10));
t122 = t114 * t112;
t110 = sin(pkin(11));
t113 = cos(pkin(11));
t117 = sin(qJ(2));
t119 = cos(qJ(2));
t121 = t119 * t110 + t117 * t113;
t120 = t117 * t110 - t119 * t113;
t118 = cos(qJ(4));
t116 = sin(qJ(4));
t115 = cos(pkin(6));
t107 = t121 * t115;
t106 = t120 * t115;
t1 = [0, t123, 0, -t111 * t106 + t114 * t121 (-t111 * t107 - t114 * t120) * t116 - t118 * t123, 0; 0, -t122, 0, t114 * t106 + t111 * t121 (t114 * t107 - t111 * t120) * t116 + t118 * t122, 0; 0, t115, 0, t120 * t112, t121 * t116 * t112 - t115 * t118, 0;];
Jg_rot  = t1;
