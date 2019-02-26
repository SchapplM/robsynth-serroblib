% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRP7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:39
% EndTime: 2019-02-26 22:12:39
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t144 = sin(pkin(6));
t147 = sin(qJ(1));
t155 = t147 * t144;
t146 = sin(qJ(2));
t154 = t147 * t146;
t148 = cos(qJ(2));
t153 = t147 * t148;
t149 = cos(qJ(1));
t152 = t149 * t144;
t151 = t149 * t146;
t150 = t149 * t148;
t145 = cos(pkin(6));
t143 = qJ(3) + pkin(11);
t142 = cos(t143);
t141 = sin(t143);
t1 = [0, t155, t145 * t153 + t151, 0 (-t145 * t154 + t150) * t141 - t142 * t155, 0; 0, -t152, -t145 * t150 + t154, 0 (t145 * t151 + t153) * t141 + t142 * t152, 0; 1, t145, -t144 * t148, 0, t144 * t146 * t141 - t145 * t142, 0;];
Jg_rot  = t1;
