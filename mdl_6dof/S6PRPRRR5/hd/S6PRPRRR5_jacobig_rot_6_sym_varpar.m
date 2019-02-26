% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:15
% EndTime: 2019-02-26 19:56:15
% DurationCPUTime: 0.03s
% Computational Cost: add. (16->9), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
t113 = sin(pkin(11));
t114 = sin(pkin(6));
t122 = t113 * t114;
t115 = cos(pkin(11));
t121 = t115 * t114;
t116 = cos(pkin(6));
t117 = sin(qJ(2));
t120 = t116 * t117;
t118 = cos(qJ(2));
t119 = t116 * t118;
t112 = qJ(4) + qJ(5);
t111 = cos(t112);
t110 = sin(t112);
t109 = t114 * t117;
t108 = -t113 * t120 + t115 * t118;
t107 = t113 * t118 + t115 * t120;
t1 = [0, t122, 0, t108, t108, t110 * t122 - (t113 * t119 + t115 * t117) * t111; 0, -t121, 0, t107, t107, -t110 * t121 - (t113 * t117 - t115 * t119) * t111; 0, t116, 0, t109, t109, t114 * t118 * t111 + t116 * t110;];
Jg_rot  = t1;
