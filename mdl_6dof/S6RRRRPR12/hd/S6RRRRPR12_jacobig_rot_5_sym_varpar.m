% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR12_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:59
% EndTime: 2019-02-26 22:36:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
t148 = sin(pkin(6));
t153 = sin(qJ(1));
t162 = t153 * t148;
t152 = sin(qJ(2));
t161 = t153 * t152;
t155 = cos(qJ(2));
t160 = t153 * t155;
t156 = cos(qJ(1));
t159 = t156 * t148;
t158 = t156 * t152;
t157 = t156 * t155;
t154 = cos(qJ(3));
t151 = sin(qJ(3));
t150 = cos(pkin(6));
t149 = cos(pkin(7));
t147 = sin(pkin(7));
t146 = -t150 * t160 - t158;
t145 = t150 * t157 - t161;
t1 = [0, t162, -t146 * t147 + t149 * t162 (-t150 * t161 + t157) * t151 + (-t146 * t149 - t147 * t162) * t154, 0, 0; 0, -t159, -t145 * t147 - t149 * t159 (t150 * t158 + t160) * t151 + (-t145 * t149 + t147 * t159) * t154, 0, 0; 1, t150, -t148 * t155 * t147 + t150 * t149, -t150 * t147 * t154 + (-t149 * t154 * t155 + t151 * t152) * t148, 0, 0;];
Jg_rot  = t1;
