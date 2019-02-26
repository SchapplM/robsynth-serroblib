% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRP1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:33
% EndTime: 2019-02-26 19:41:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
t138 = sin(pkin(11));
t144 = cos(pkin(6));
t156 = t138 * t144;
t139 = sin(pkin(7));
t140 = sin(pkin(6));
t155 = t139 * t140;
t154 = t139 * t144;
t143 = cos(pkin(7));
t153 = t140 * t143;
t141 = cos(pkin(12));
t152 = t141 * t143;
t142 = cos(pkin(11));
t151 = t142 * t144;
t137 = sin(pkin(12));
t133 = -t138 * t137 + t141 * t151;
t150 = -t133 * t143 + t142 * t155;
t135 = -t142 * t137 - t141 * t156;
t149 = t135 * t143 + t138 * t155;
t148 = cos(qJ(3));
t147 = cos(qJ(4));
t146 = sin(qJ(3));
t145 = sin(qJ(4));
t136 = -t137 * t156 + t142 * t141;
t134 = t137 * t151 + t138 * t141;
t132 = -t141 * t155 + t144 * t143;
t131 = -t135 * t139 + t138 * t153;
t130 = -t133 * t139 - t142 * t153;
t1 = [0, 0, t131, t136 * t146 - t149 * t148 (t136 * t148 + t149 * t146) * t145 - t131 * t147, 0; 0, 0, t130, t134 * t146 + t150 * t148 (t134 * t148 - t150 * t146) * t145 - t130 * t147, 0; 0, 0, t132, -t148 * t154 + (t137 * t146 - t148 * t152) * t140 (t146 * t154 + (t137 * t148 + t146 * t152) * t140) * t145 - t132 * t147, 0;];
Jg_rot  = t1;
