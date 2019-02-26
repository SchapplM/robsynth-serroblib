% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:07
% EndTime: 2019-02-26 22:34:07
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t135 = sin(pkin(6));
t139 = cos(qJ(2));
t147 = t135 * t139;
t138 = sin(qJ(1));
t146 = t138 * t135;
t137 = sin(qJ(2));
t145 = t138 * t137;
t144 = t138 * t139;
t140 = cos(qJ(1));
t143 = t140 * t135;
t142 = t140 * t137;
t141 = t140 * t139;
t136 = cos(pkin(6));
t134 = qJ(3) + qJ(4) + pkin(12);
t133 = cos(t134);
t132 = sin(t134);
t131 = t136 * t144 + t142;
t130 = -t136 * t141 + t145;
t1 = [0, t146, t131, t131, 0 (-t136 * t145 + t141) * t132 - t133 * t146; 0, -t143, t130, t130, 0 (t136 * t142 + t144) * t132 + t133 * t143; 1, t136, -t147, -t147, 0, t135 * t137 * t132 - t136 * t133;];
Jg_rot  = t1;
