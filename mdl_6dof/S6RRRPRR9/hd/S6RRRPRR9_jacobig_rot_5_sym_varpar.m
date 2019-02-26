% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR9_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:28
% EndTime: 2019-02-26 22:20:28
% DurationCPUTime: 0.07s
% Computational Cost: add. (25->16), mult. (72->33), div. (0->0), fcn. (105->12), ass. (0->25)
t151 = sin(pkin(6));
t157 = sin(qJ(1));
t167 = t157 * t151;
t156 = sin(qJ(2));
t166 = t157 * t156;
t159 = cos(qJ(2));
t165 = t157 * t159;
t160 = cos(qJ(1));
t164 = t160 * t151;
t163 = t160 * t156;
t162 = t160 * t159;
t149 = sin(pkin(13));
t152 = cos(pkin(13));
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t161 = -t149 * t155 + t152 * t158;
t154 = cos(pkin(6));
t153 = cos(pkin(7));
t150 = sin(pkin(7));
t148 = -t149 * t158 - t152 * t155;
t147 = -t154 * t165 - t163;
t146 = t154 * t162 - t166;
t145 = t161 * t153;
t144 = t161 * t150;
t1 = [0, t167, -t147 * t150 + t153 * t167, 0 -(-t154 * t166 + t162) * t148 - t147 * t145 - t144 * t167, 0; 0, -t164, -t146 * t150 - t153 * t164, 0 -(t154 * t163 + t165) * t148 - t146 * t145 + t144 * t164, 0; 1, t154, -t150 * t151 * t159 + t153 * t154, 0, -t154 * t144 + (-t145 * t159 - t148 * t156) * t151, 0;];
Jg_rot  = t1;
