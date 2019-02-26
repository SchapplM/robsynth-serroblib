% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:42
% EndTime: 2019-02-26 20:04:42
% DurationCPUTime: 0.03s
% Computational Cost: add. (26->10), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->17)
t114 = sin(pkin(11));
t115 = sin(pkin(6));
t123 = t114 * t115;
t116 = cos(pkin(11));
t122 = t116 * t115;
t117 = cos(pkin(6));
t118 = sin(qJ(2));
t121 = t117 * t118;
t119 = cos(qJ(2));
t120 = t117 * t119;
t113 = qJ(3) + pkin(12);
t112 = cos(t113);
t111 = sin(t113);
t110 = t115 * t118 * t111 - t117 * t112;
t109 = (-t114 * t121 + t116 * t119) * t111 - t112 * t123;
t108 = (t114 * t119 + t116 * t121) * t111 + t112 * t122;
t1 = [0, t123, t114 * t120 + t116 * t118, 0, t109, t109; 0, -t122, t114 * t118 - t116 * t120, 0, t108, t108; 0, t117, -t115 * t119, 0, t110, t110;];
Jg_rot  = t1;
