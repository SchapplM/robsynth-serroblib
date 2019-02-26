% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR3_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.04s
% Computational Cost: add. (25->16), mult. (72->35), div. (0->0), fcn. (105->12), ass. (0->23)
t133 = sin(pkin(12));
t135 = sin(pkin(6));
t148 = t133 * t135;
t137 = cos(pkin(12));
t147 = t137 * t135;
t139 = cos(pkin(6));
t141 = sin(qJ(2));
t146 = t139 * t141;
t143 = cos(qJ(2));
t145 = t139 * t143;
t132 = sin(pkin(13));
t136 = cos(pkin(13));
t140 = sin(qJ(3));
t142 = cos(qJ(3));
t144 = -t132 * t140 + t136 * t142;
t138 = cos(pkin(7));
t134 = sin(pkin(7));
t131 = -t142 * t132 - t140 * t136;
t130 = -t133 * t145 - t137 * t141;
t129 = -t133 * t141 + t137 * t145;
t128 = t144 * t138;
t127 = t144 * t134;
t1 = [0, t148, -t130 * t134 + t138 * t148, 0 -(-t133 * t146 + t137 * t143) * t131 - t130 * t128 - t127 * t148, 0; 0, -t147, -t129 * t134 - t138 * t147, 0 -(t133 * t143 + t137 * t146) * t131 - t129 * t128 + t127 * t147, 0; 0, t139, -t135 * t143 * t134 + t139 * t138, 0, -t139 * t127 + (-t128 * t143 - t131 * t141) * t135, 0;];
Jg_rot  = t1;
