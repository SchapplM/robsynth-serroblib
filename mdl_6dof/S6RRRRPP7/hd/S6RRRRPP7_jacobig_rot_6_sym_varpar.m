% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:55
% EndTime: 2019-02-26 22:28:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t141 = sin(pkin(6));
t145 = sin(qJ(1));
t154 = t145 * t141;
t144 = sin(qJ(2));
t153 = t145 * t144;
t147 = cos(qJ(2));
t152 = t145 * t147;
t148 = cos(qJ(1));
t151 = t148 * t141;
t150 = t148 * t144;
t149 = t148 * t147;
t146 = cos(qJ(3));
t143 = sin(qJ(3));
t142 = cos(pkin(6));
t1 = [0, t154, t142 * t152 + t150 (-t142 * t153 + t149) * t143 - t146 * t154, 0, 0; 0, -t151, -t142 * t149 + t153 (t142 * t150 + t152) * t143 + t146 * t151, 0, 0; 1, t142, -t141 * t147, t141 * t144 * t143 - t142 * t146, 0, 0;];
Jg_rot  = t1;
