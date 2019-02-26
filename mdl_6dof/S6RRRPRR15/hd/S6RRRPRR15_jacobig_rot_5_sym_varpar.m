% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR15_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:25
% EndTime: 2019-02-26 22:24:25
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
t140 = sin(pkin(6));
t145 = sin(qJ(1));
t154 = t145 * t140;
t144 = sin(qJ(2));
t153 = t145 * t144;
t147 = cos(qJ(2));
t152 = t145 * t147;
t148 = cos(qJ(1));
t151 = t148 * t140;
t150 = t148 * t144;
t149 = t148 * t147;
t146 = cos(qJ(3));
t143 = sin(qJ(3));
t142 = cos(pkin(6));
t141 = cos(pkin(7));
t139 = sin(pkin(7));
t138 = -t142 * t152 - t150;
t137 = t142 * t149 - t153;
t1 = [0, t154, -t138 * t139 + t141 * t154, 0 (-t142 * t153 + t149) * t146 + (t138 * t141 + t139 * t154) * t143, 0; 0, -t151, -t137 * t139 - t141 * t151, 0 (t142 * t150 + t152) * t146 + (t137 * t141 - t139 * t151) * t143, 0; 1, t142, -t140 * t147 * t139 + t142 * t141, 0, t142 * t139 * t143 + (t141 * t143 * t147 + t144 * t146) * t140, 0;];
Jg_rot  = t1;
