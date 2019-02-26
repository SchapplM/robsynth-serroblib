% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRPR2_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
t109 = sin(pkin(11));
t115 = cos(pkin(6));
t121 = t109 * t115;
t110 = sin(pkin(7));
t111 = sin(pkin(6));
t120 = t110 * t111;
t114 = cos(pkin(7));
t119 = t111 * t114;
t113 = cos(pkin(11));
t118 = t113 * t115;
t117 = cos(qJ(3));
t116 = sin(qJ(3));
t112 = cos(pkin(12));
t108 = sin(pkin(12));
t107 = -t113 * t108 - t112 * t121;
t106 = -t109 * t108 + t112 * t118;
t1 = [0, 0, -t107 * t110 + t109 * t119 (-t108 * t121 + t113 * t112) * t116 + (-t107 * t114 - t109 * t120) * t117, 0, 0; 0, 0, -t106 * t110 - t113 * t119 (t108 * t118 + t109 * t112) * t116 + (-t106 * t114 + t113 * t120) * t117, 0, 0; 0, 0, -t112 * t120 + t115 * t114, -t115 * t110 * t117 + (-t112 * t114 * t117 + t108 * t116) * t111, 0, 0;];
Jg_rot  = t1;
