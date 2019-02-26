% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR8_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:42
% EndTime: 2019-02-26 20:14:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
t140 = sin(pkin(12));
t142 = sin(pkin(6));
t154 = t140 * t142;
t141 = sin(pkin(7));
t153 = t141 * t142;
t143 = cos(pkin(12));
t152 = t143 * t142;
t145 = cos(pkin(6));
t147 = sin(qJ(2));
t151 = t145 * t147;
t149 = cos(qJ(2));
t150 = t145 * t149;
t148 = cos(qJ(3));
t146 = sin(qJ(3));
t144 = cos(pkin(7));
t139 = -t140 * t150 - t143 * t147;
t138 = -t140 * t147 + t143 * t150;
t1 = [0, t154, -t139 * t141 + t144 * t154 (-t140 * t151 + t143 * t149) * t146 + (-t139 * t144 - t140 * t153) * t148, 0, 0; 0, -t152, -t138 * t141 - t144 * t152 (t140 * t149 + t143 * t151) * t146 + (-t138 * t144 + t141 * t152) * t148, 0, 0; 0, t145, t145 * t144 - t149 * t153, -t145 * t141 * t148 + (-t144 * t148 * t149 + t146 * t147) * t142, 0, 0;];
Jg_rot  = t1;
