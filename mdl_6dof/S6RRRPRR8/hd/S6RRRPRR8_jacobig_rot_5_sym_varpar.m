% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR8_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:48
% EndTime: 2019-02-26 22:19:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t112 = sin(pkin(6));
t115 = sin(qJ(1));
t123 = t115 * t112;
t114 = sin(qJ(2));
t122 = t115 * t114;
t116 = cos(qJ(2));
t121 = t115 * t116;
t117 = cos(qJ(1));
t120 = t117 * t112;
t119 = t117 * t114;
t118 = t117 * t116;
t113 = cos(pkin(6));
t111 = qJ(3) + pkin(12);
t110 = cos(t111);
t109 = sin(t111);
t1 = [0, t123, t113 * t121 + t119, 0 (-t113 * t122 + t118) * t109 - t110 * t123, 0; 0, -t120, -t113 * t118 + t122, 0 (t113 * t119 + t121) * t109 + t110 * t120, 0; 1, t113, -t112 * t116, 0, t112 * t114 * t109 - t113 * t110, 0;];
Jg_rot  = t1;
