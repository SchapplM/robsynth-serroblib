% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRP6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:03
% EndTime: 2019-02-26 22:12:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t114 = sin(pkin(6));
t117 = sin(qJ(1));
t125 = t117 * t114;
t116 = sin(qJ(2));
t124 = t117 * t116;
t118 = cos(qJ(2));
t123 = t117 * t118;
t119 = cos(qJ(1));
t122 = t119 * t114;
t121 = t119 * t116;
t120 = t119 * t118;
t115 = cos(pkin(6));
t113 = qJ(3) + pkin(11);
t112 = cos(t113);
t111 = sin(t113);
t1 = [0, t125, t115 * t123 + t121, 0 (-t115 * t124 + t120) * t111 - t112 * t125, 0; 0, -t122, -t115 * t120 + t124, 0 (t115 * t121 + t123) * t111 + t112 * t122, 0; 1, t115, -t114 * t118, 0, t114 * t116 * t111 - t115 * t112, 0;];
Jg_rot  = t1;
