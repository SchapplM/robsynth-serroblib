% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR10_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:18
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t111 = sin(pkin(6));
t114 = sin(qJ(1));
t122 = t114 * t111;
t113 = sin(qJ(2));
t121 = t114 * t113;
t115 = cos(qJ(2));
t120 = t114 * t115;
t116 = cos(qJ(1));
t119 = t116 * t111;
t118 = t116 * t113;
t117 = t116 * t115;
t112 = cos(pkin(6));
t110 = pkin(12) + qJ(4);
t109 = cos(t110);
t108 = sin(t110);
t1 = [0, t122, 0, t112 * t120 + t118 (-t112 * t121 + t117) * t108 - t109 * t122, 0; 0, -t119, 0, -t112 * t117 + t121 (t112 * t118 + t120) * t108 + t109 * t119, 0; 1, t112, 0, -t111 * t115, t111 * t113 * t108 - t112 * t109, 0;];
Jg_rot  = t1;
