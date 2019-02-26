% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:43
% DurationCPUTime: 0.06s
% Computational Cost: add. (31->11), mult. (70->24), div. (0->0), fcn. (106->10), ass. (0->21)
t143 = sin(pkin(11));
t144 = sin(pkin(6));
t153 = t143 * t144;
t146 = cos(pkin(11));
t152 = t146 * t144;
t142 = sin(pkin(12));
t145 = cos(pkin(12));
t148 = sin(qJ(2));
t149 = cos(qJ(2));
t151 = t149 * t142 + t148 * t145;
t150 = t148 * t142 - t149 * t145;
t147 = cos(pkin(6));
t141 = qJ(4) + qJ(5);
t140 = cos(t141);
t139 = sin(t141);
t136 = t151 * t147;
t135 = t150 * t147;
t134 = t150 * t144;
t133 = -t143 * t135 + t146 * t151;
t132 = t146 * t135 + t143 * t151;
t1 = [0, t153, 0, t133, t133 (-t143 * t136 - t146 * t150) * t139 - t140 * t153; 0, -t152, 0, t132, t132 (t146 * t136 - t143 * t150) * t139 + t140 * t152; 0, t147, 0, t134, t134, t151 * t139 * t144 - t147 * t140;];
Jg_rot  = t1;
