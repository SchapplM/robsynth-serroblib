% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPPR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:50
% EndTime: 2019-02-26 22:05:50
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t113 = sin(pkin(6));
t116 = sin(qJ(1));
t124 = t116 * t113;
t115 = sin(qJ(2));
t123 = t116 * t115;
t117 = cos(qJ(2));
t122 = t116 * t117;
t118 = cos(qJ(1));
t121 = t118 * t113;
t120 = t118 * t115;
t119 = t118 * t117;
t114 = cos(pkin(6));
t112 = qJ(3) + pkin(11);
t111 = cos(t112);
t110 = sin(t112);
t1 = [0, t124, t114 * t122 + t120, 0, 0 (-t114 * t123 + t119) * t110 - t111 * t124; 0, -t121, -t114 * t119 + t123, 0, 0 (t114 * t120 + t122) * t110 + t111 * t121; 1, t114, -t113 * t117, 0, 0, t113 * t115 * t110 - t114 * t111;];
Jg_rot  = t1;
