% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:23
% EndTime: 2019-02-26 21:42:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t114 = sin(qJ(2));
t117 = cos(qJ(1));
t123 = t114 * t117;
t112 = sin(pkin(6));
t115 = sin(qJ(1));
t122 = t115 * t112;
t121 = t115 * t114;
t116 = cos(qJ(2));
t120 = t115 * t116;
t119 = t117 * t112;
t118 = t117 * t116;
t113 = cos(pkin(6));
t111 = pkin(11) + qJ(4);
t110 = cos(t111);
t109 = sin(t111);
t1 = [0, t122, 0, t113 * t120 + t123, 0 (-t113 * t121 + t118) * t109 - t110 * t122; 0, -t119, 0, -t113 * t118 + t121, 0 (t113 * t123 + t120) * t109 + t110 * t119; 1, t113, 0, -t112 * t116, 0, t109 * t112 * t114 - t113 * t110;];
Jg_rot  = t1;
