% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:54
% EndTime: 2019-02-26 20:16:54
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->9), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
t137 = sin(pkin(11));
t138 = sin(pkin(6));
t148 = t137 * t138;
t139 = cos(pkin(11));
t147 = t139 * t138;
t140 = cos(pkin(6));
t142 = sin(qJ(2));
t146 = t140 * t142;
t144 = cos(qJ(2));
t145 = t140 * t144;
t143 = cos(qJ(3));
t141 = sin(qJ(3));
t136 = t138 * t142 * t141 - t140 * t143;
t135 = (-t137 * t146 + t139 * t144) * t141 - t143 * t148;
t134 = (t137 * t144 + t139 * t146) * t141 + t143 * t147;
t1 = [0, t148, t137 * t145 + t139 * t142, t135, t135, 0; 0, -t147, t137 * t142 - t139 * t145, t134, t134, 0; 0, t140, -t138 * t144, t136, t136, 0;];
Jg_rot  = t1;
