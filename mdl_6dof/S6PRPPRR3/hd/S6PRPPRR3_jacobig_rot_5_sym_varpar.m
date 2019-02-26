% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPPRR3_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (22->19), div. (0->0), fcn. (35->8), ass. (0->11)
t79 = cos(pkin(6));
t80 = sin(qJ(2));
t83 = t79 * t80;
t81 = cos(qJ(2));
t82 = t79 * t81;
t78 = cos(pkin(10));
t77 = cos(pkin(11));
t76 = sin(pkin(6));
t75 = sin(pkin(10));
t74 = sin(pkin(11));
t1 = [0, t75 * t76, 0, 0 (-t75 * t83 + t78 * t81) * t74 - (t75 * t82 + t78 * t80) * t77, 0; 0, -t78 * t76, 0, 0 (t75 * t81 + t78 * t83) * t74 - (t75 * t80 - t78 * t82) * t77, 0; 0, t79, 0, 0 (t74 * t80 + t77 * t81) * t76, 0;];
Jg_rot  = t1;
