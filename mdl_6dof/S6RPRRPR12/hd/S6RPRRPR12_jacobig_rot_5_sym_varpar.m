% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR12_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:16
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
t141 = sin(pkin(6));
t146 = sin(qJ(1));
t154 = t141 * t146;
t148 = cos(qJ(1));
t153 = t141 * t148;
t139 = sin(pkin(12));
t152 = t146 * t139;
t142 = cos(pkin(12));
t151 = t146 * t142;
t150 = t148 * t139;
t149 = t148 * t142;
t147 = cos(qJ(3));
t145 = sin(qJ(3));
t144 = cos(pkin(6));
t143 = cos(pkin(7));
t140 = sin(pkin(7));
t138 = -t144 * t151 - t150;
t137 = t144 * t149 - t152;
t1 = [0, 0, -t138 * t140 + t143 * t154 (-t144 * t152 + t149) * t145 + (-t138 * t143 - t140 * t154) * t147, 0, 0; 0, 0, -t137 * t140 - t143 * t153 (t144 * t150 + t151) * t145 + (-t137 * t143 + t140 * t153) * t147, 0, 0; 1, 0, -t141 * t142 * t140 + t144 * t143, -t144 * t140 * t147 + (-t142 * t143 * t147 + t139 * t145) * t141, 0, 0;];
Jg_rot  = t1;
