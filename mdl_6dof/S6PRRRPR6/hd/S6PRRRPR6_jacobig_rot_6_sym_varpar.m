% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:35
% EndTime: 2019-02-26 20:13:35
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->12), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
t114 = sin(pkin(11));
t115 = sin(pkin(6));
t125 = t114 * t115;
t116 = cos(pkin(11));
t124 = t116 * t115;
t117 = cos(pkin(6));
t119 = sin(qJ(2));
t123 = t117 * t119;
t121 = cos(qJ(2));
t122 = t117 * t121;
t120 = cos(qJ(3));
t118 = sin(qJ(3));
t113 = t115 * t119 * t118 - t117 * t120;
t112 = (-t114 * t123 + t116 * t121) * t118 - t120 * t125;
t111 = (t114 * t121 + t116 * t123) * t118 + t120 * t124;
t1 = [0, t125, t114 * t122 + t116 * t119, t112, 0, -t112; 0, -t124, t114 * t119 - t116 * t122, t111, 0, -t111; 0, t117, -t115 * t121, t113, 0, -t113;];
Jg_rot  = t1;
