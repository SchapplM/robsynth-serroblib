% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:27
% EndTime: 2019-02-26 20:46:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t76 = sin(qJ(5));
t77 = sin(qJ(1));
t83 = t77 * t76;
t78 = cos(qJ(5));
t82 = t77 * t78;
t79 = cos(qJ(1));
t81 = t79 * t76;
t80 = t79 * t78;
t75 = pkin(9) + qJ(3);
t74 = cos(t75);
t73 = sin(t75);
t72 = -t73 * t83 + t80;
t71 = t73 * t82 + t81;
t70 = t73 * t81 + t82;
t69 = t73 * t80 - t83;
t1 = [t72, 0, t74 * t81, 0, t69, 0; t70, 0, t74 * t83, 0, t71, 0; 0, 0, t73 * t76, 0, -t74 * t78, 0; -t71, 0, t74 * t80, 0, -t70, 0; t69, 0, t74 * t82, 0, t72, 0; 0, 0, t73 * t78, 0, t74 * t76, 0; -t77 * t74, 0, -t79 * t73, 0, 0, 0; t79 * t74, 0, -t77 * t73, 0, 0, 0; 0, 0, t74, 0, 0, 0;];
JR_rot  = t1;
