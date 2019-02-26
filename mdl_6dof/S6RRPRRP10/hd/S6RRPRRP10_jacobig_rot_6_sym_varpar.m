% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:04
% EndTime: 2019-02-26 21:51:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t143 = sin(pkin(6));
t146 = sin(qJ(1));
t154 = t146 * t143;
t145 = sin(qJ(2));
t153 = t146 * t145;
t147 = cos(qJ(2));
t152 = t146 * t147;
t148 = cos(qJ(1));
t151 = t148 * t143;
t150 = t148 * t145;
t149 = t148 * t147;
t144 = cos(pkin(6));
t142 = pkin(11) + qJ(4);
t141 = cos(t142);
t140 = sin(t142);
t1 = [0, t154, 0, t144 * t152 + t150 (-t144 * t153 + t149) * t140 - t141 * t154, 0; 0, -t151, 0, -t144 * t149 + t153 (t144 * t150 + t152) * t140 + t141 * t151, 0; 1, t144, 0, -t143 * t147, t143 * t145 * t140 - t144 * t141, 0;];
Jg_rot  = t1;
