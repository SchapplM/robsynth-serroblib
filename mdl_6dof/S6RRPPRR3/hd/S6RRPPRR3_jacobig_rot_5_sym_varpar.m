% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Jg_rot = S6RRPPRR3_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->6), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
t101 = cos(pkin(11));
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t99 = sin(pkin(11));
t107 = t101 * t105 - t103 * t99;
t106 = cos(qJ(1));
t104 = sin(qJ(1));
t102 = cos(pkin(6));
t100 = sin(pkin(6));
t98 = -t103 * t101 - t105 * t99;
t97 = t107 * t102;
t1 = [0, t104 * t100, 0, 0, t104 * t97 - t106 * t98, 0; 0, -t106 * t100, 0, 0, -t104 * t98 - t106 * t97, 0; 1, t102, 0, 0, -t107 * t100, 0;];
Jg_rot  = t1;
