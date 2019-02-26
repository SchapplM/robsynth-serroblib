% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:01
% EndTime: 2019-02-26 21:43:01
% DurationCPUTime: 0.08s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t107 = sin(pkin(6));
t110 = sin(qJ(1));
t118 = t110 * t107;
t109 = sin(qJ(2));
t117 = t110 * t109;
t111 = cos(qJ(2));
t116 = t110 * t111;
t112 = cos(qJ(1));
t115 = t112 * t107;
t114 = t112 * t109;
t113 = t112 * t111;
t108 = cos(pkin(6));
t106 = pkin(11) + qJ(4);
t105 = cos(t106);
t104 = sin(t106);
t1 = [0, t118, 0, t108 * t116 + t114, 0 (-t108 * t117 + t113) * t105 + t104 * t118; 0, -t115, 0, -t108 * t113 + t117, 0 (t108 * t114 + t116) * t105 - t104 * t115; 1, t108, 0, -t107 * t111, 0, t107 * t109 * t105 + t108 * t104;];
Jg_rot  = t1;
