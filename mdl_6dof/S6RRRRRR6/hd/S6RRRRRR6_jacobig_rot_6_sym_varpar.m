% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:04
% EndTime: 2019-02-26 22:50:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->11), mult. (46->18), div. (0->0), fcn. (75->8), ass. (0->22)
t147 = sin(pkin(6));
t151 = cos(qJ(2));
t159 = t147 * t151;
t150 = sin(qJ(1));
t158 = t150 * t147;
t149 = sin(qJ(2));
t157 = t150 * t149;
t156 = t150 * t151;
t152 = cos(qJ(1));
t155 = t152 * t147;
t154 = t152 * t149;
t153 = t152 * t151;
t148 = cos(pkin(6));
t146 = qJ(3) + qJ(4);
t145 = cos(t146);
t144 = sin(t146);
t143 = t148 * t156 + t154;
t142 = -t148 * t153 + t157;
t141 = t147 * t149 * t144 - t148 * t145;
t140 = (-t148 * t157 + t153) * t144 - t145 * t158;
t139 = (t148 * t154 + t156) * t144 + t145 * t155;
t1 = [0, t158, t143, t143, t140, t140; 0, -t155, t142, t142, t139, t139; 1, t148, -t159, -t159, t141, t141;];
Jg_rot  = t1;