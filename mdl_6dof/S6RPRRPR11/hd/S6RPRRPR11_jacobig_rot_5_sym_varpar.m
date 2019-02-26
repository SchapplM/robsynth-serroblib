% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR11_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:40
% EndTime: 2019-02-26 21:06:40
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
t151 = sin(pkin(6));
t156 = sin(qJ(1));
t164 = t151 * t156;
t158 = cos(qJ(1));
t163 = t151 * t158;
t149 = sin(pkin(12));
t162 = t156 * t149;
t152 = cos(pkin(12));
t161 = t156 * t152;
t160 = t158 * t149;
t159 = t158 * t152;
t157 = cos(qJ(3));
t155 = sin(qJ(3));
t154 = cos(pkin(6));
t153 = cos(pkin(7));
t150 = sin(pkin(7));
t148 = -t154 * t161 - t160;
t147 = t154 * t159 - t162;
t1 = [0, 0, -t148 * t150 + t153 * t164 (-t154 * t162 + t159) * t155 + (-t148 * t153 - t150 * t164) * t157, 0, 0; 0, 0, -t147 * t150 - t153 * t163 (t154 * t160 + t161) * t155 + (-t147 * t153 + t150 * t163) * t157, 0, 0; 1, 0, -t151 * t152 * t150 + t154 * t153, -t154 * t150 * t157 + (-t152 * t153 * t157 + t149 * t155) * t151, 0, 0;];
Jg_rot  = t1;
