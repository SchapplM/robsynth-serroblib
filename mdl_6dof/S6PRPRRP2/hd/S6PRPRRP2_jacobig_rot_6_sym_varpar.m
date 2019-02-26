% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:04
% EndTime: 2019-02-26 19:51:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t136 = sin(pkin(10));
t137 = sin(pkin(6));
t148 = t136 * t137;
t139 = cos(pkin(10));
t147 = t139 * t137;
t135 = sin(pkin(11));
t138 = cos(pkin(11));
t142 = sin(qJ(2));
t144 = cos(qJ(2));
t146 = t144 * t135 + t142 * t138;
t145 = t142 * t135 - t144 * t138;
t143 = cos(qJ(4));
t141 = sin(qJ(4));
t140 = cos(pkin(6));
t132 = t146 * t140;
t131 = t145 * t140;
t1 = [0, t148, 0, -t136 * t131 + t139 * t146 (-t136 * t132 - t139 * t145) * t141 - t143 * t148, 0; 0, -t147, 0, t139 * t131 + t136 * t146 (t139 * t132 - t136 * t145) * t141 + t143 * t147, 0; 0, t140, 0, t145 * t137, t146 * t141 * t137 - t140 * t143, 0;];
Jg_rot  = t1;
