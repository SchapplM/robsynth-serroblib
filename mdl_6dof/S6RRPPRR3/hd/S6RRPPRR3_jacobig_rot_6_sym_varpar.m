% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPPRR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:34
% EndTime: 2019-02-26 21:29:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->11), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->18)
t139 = sin(pkin(6));
t143 = sin(qJ(1));
t149 = t143 * t139;
t145 = cos(qJ(1));
t148 = t145 * t139;
t138 = sin(pkin(11));
t140 = cos(pkin(11));
t142 = sin(qJ(2));
t144 = cos(qJ(2));
t147 = t144 * t138 + t142 * t140;
t146 = t142 * t138 - t144 * t140;
t141 = cos(pkin(6));
t137 = pkin(12) + qJ(5);
t136 = cos(t137);
t135 = sin(t137);
t132 = t147 * t141;
t131 = t146 * t141;
t1 = [0, t149, 0, 0, -t143 * t131 + t145 * t147 (-t143 * t132 - t145 * t146) * t135 - t136 * t149; 0, -t148, 0, 0, t145 * t131 + t143 * t147 (t145 * t132 - t143 * t146) * t135 + t136 * t148; 1, t141, 0, 0, t146 * t139, t147 * t135 * t139 - t141 * t136;];
Jg_rot  = t1;
