% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:19
% EndTime: 2019-02-26 21:56:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (28->10), mult. (78->24), div. (0->0), fcn. (117->10), ass. (0->20)
t146 = sin(pkin(6));
t151 = sin(qJ(1));
t158 = t151 * t146;
t154 = cos(qJ(1));
t157 = t154 * t146;
t145 = sin(pkin(12));
t147 = cos(pkin(12));
t150 = sin(qJ(2));
t153 = cos(qJ(2));
t156 = t153 * t145 + t150 * t147;
t155 = t150 * t145 - t153 * t147;
t152 = cos(qJ(4));
t149 = sin(qJ(4));
t148 = cos(pkin(6));
t142 = t156 * t148;
t141 = t155 * t148;
t140 = t156 * t149 * t146 - t148 * t152;
t139 = (-t151 * t142 - t154 * t155) * t149 - t152 * t158;
t138 = (t154 * t142 - t151 * t155) * t149 + t152 * t157;
t1 = [0, t158, 0, -t151 * t141 + t154 * t156, t139, t139; 0, -t157, 0, t154 * t141 + t151 * t156, t138, t138; 1, t148, 0, t155 * t146, t140, t140;];
Jg_rot  = t1;
