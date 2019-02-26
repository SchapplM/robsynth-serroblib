% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR12_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:08
% EndTime: 2019-02-26 21:21:08
% DurationCPUTime: 0.05s
% Computational Cost: add. (23->16), mult. (67->37), div. (0->0), fcn. (96->12), ass. (0->24)
t142 = sin(pkin(6));
t148 = sin(qJ(1));
t156 = t142 * t148;
t150 = cos(qJ(1));
t155 = t142 * t150;
t139 = sin(pkin(14));
t154 = t148 * t139;
t143 = cos(pkin(14));
t153 = t148 * t143;
t152 = t150 * t139;
t151 = t150 * t143;
t149 = cos(qJ(3));
t147 = sin(qJ(3));
t146 = cos(pkin(6));
t145 = cos(pkin(7));
t144 = cos(pkin(8));
t141 = sin(pkin(7));
t140 = sin(pkin(8));
t138 = -t146 * t153 - t152;
t137 = t146 * t151 - t154;
t136 = -t142 * t143 * t141 + t146 * t145;
t135 = -t138 * t141 + t145 * t156;
t134 = -t137 * t141 - t145 * t155;
t1 = [0, 0, t135 -(-(-t146 * t154 + t151) * t147 + (t138 * t145 + t141 * t156) * t149) * t140 + t135 * t144, 0, 0; 0, 0, t134 -(-(t146 * t152 + t153) * t147 + (t137 * t145 - t141 * t155) * t149) * t140 + t134 * t144, 0, 0; 1, 0, t136 -(t146 * t141 * t149 + (t143 * t145 * t149 - t139 * t147) * t142) * t140 + t136 * t144, 0, 0;];
Jg_rot  = t1;
