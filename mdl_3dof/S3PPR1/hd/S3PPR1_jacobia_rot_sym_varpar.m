% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S3PPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% Ja_rot [3x3]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S3PPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_jacobia_rot_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3PPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_jacobia_rot_sym_varpar: pkin has to be [3x1] (double)');
%% Function calls
if link_index == 0
	Ja_rot=S3PPR1_jacobia_rot_0_sym_varpar(qJ, pkin);
elseif link_index == 1
	Ja_rot=S3PPR1_jacobia_rot_1_sym_varpar(qJ, pkin);
elseif link_index == 2
	Ja_rot=S3PPR1_jacobia_rot_2_sym_varpar(qJ, pkin);
elseif link_index == 3
	Ja_rot=S3PPR1_jacobia_rot_3_sym_varpar(qJ, pkin);
else
	Ja_rot=NaN(3,3);
end