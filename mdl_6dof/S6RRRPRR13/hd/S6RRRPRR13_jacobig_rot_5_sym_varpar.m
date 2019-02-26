% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR13_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
t147 = sin(pkin(6));
t152 = sin(qJ(1));
t161 = t152 * t147;
t151 = sin(qJ(2));
t160 = t152 * t151;
t154 = cos(qJ(2));
t159 = t152 * t154;
t155 = cos(qJ(1));
t158 = t155 * t147;
t157 = t155 * t151;
t156 = t155 * t154;
t153 = cos(qJ(3));
t150 = sin(qJ(3));
t149 = cos(pkin(6));
t148 = cos(pkin(7));
t146 = sin(pkin(7));
t145 = -t149 * t159 - t157;
t144 = t149 * t156 - t160;
t1 = [0, t161, -t145 * t146 + t148 * t161, 0 (-t149 * t160 + t156) * t150 + (-t145 * t148 - t146 * t161) * t153, 0; 0, -t158, -t144 * t146 - t148 * t158, 0 (t149 * t157 + t159) * t150 + (-t144 * t148 + t146 * t158) * t153, 0; 1, t149, -t147 * t154 * t146 + t149 * t148, 0, -t149 * t146 * t153 + (-t148 * t153 * t154 + t150 * t151) * t147, 0;];
Jg_rot  = t1;
